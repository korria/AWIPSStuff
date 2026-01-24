#!/awips2/python/bin/python
"""
alertMon.py -  This progam monitors a serial line for ALERT GAUGE data.
Single-port version - run multiple instances for multiple ports as a daemon

Usage:
    ./alertMon.py -p 1    # Monitor /dev/ttyUSB0
    ./alertMon.py -p 2    # Monitor /dev/ttyUSB1
    ./alertMon.py -p 3    # Monitor /dev/ttyUSB2
    ./alertMon.py -p 4    # Monitor /dev/ttyUSB3
"""

import argparse
import os
import signal
import socket
import shutil
import syslog
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, '/ldad/localapps/alert/lib')

try:
    import serial
except ImportError:
    sys.exit(1)

##########################################################################
# CONFIGURATION
##########################################################################
ALERTHOME = Path("/ldad/localapps/alert")
ALERTLOGDIR = Path("/data/logs/ldad")
CONFIG_FILE = ALERTHOME / "config" / "AlertStationInfo.txt"

PILID = "BOIRRXBOI"
WFOSITE = "BOISE ID"
RR5DATAFILE = ALERTHOME / "data" / "rr5data.txt"

DEFAULT_BAUD = 9600
PORT_MAP = {
    1: "/dev/ttyUSB0",
    2: "/dev/ttyUSB1",
    3: "/dev/ttyUSB2",
    4: "/dev/ttyUSB3",
}

# ALERT format masks
ADFMASK1 = 0x0F
BDFMASK1 = 0x3F
BDFMASK2 = 0x3E
BDFMASK3 = 0x01
EIFMASK1 = 0x3F
EIFMASK2 = 0x80
EIFMASK3 = 0x7F
EIFMASK5 = 0x03

# Global state
STOP = False
logfd = None
DEBUG = False
USE_SYSLOG = False
PGMNAME = "alertMon"


def log_message(message, level=syslog.LOG_INFO):
    """Log to file and optionally syslog."""
    global logfd
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    if logfd:
        logfd.write("{}: {}\n".format(timestamp, message))
        logfd.flush()
    
    if USE_SYSLOG:
        syslog.syslog(level, message)
    
    if DEBUG:
        print("{}: {}".format(timestamp, message))


def open_log_file(logdir, pgmname):
    """Open a new log file."""
    global logfd
    
    now = datetime.now()
    log_subdir = logdir / now.strftime('%Y%m%d')
    log_subdir.mkdir(parents=True, exist_ok=True)
    
    filename = "{}_{}-{}.log".format(pgmname, os.getpid(), now.strftime('%H%M%S'))
    logpath = log_subdir / filename
    
    logfd = open(logpath, 'a')
    return logpath


def rotate_log(signum, frame):
    """Signal handler for SIGHUP - rotate log file."""
    global logfd
    log_message("Received SIGHUP - rotating log file")
    if logfd:
        logfd.close()
    open_log_file(ALERTLOGDIR, PGMNAME)


def exit_handler(signum, frame):
    """Signal handler for clean exit."""
    global STOP
    log_message("Received exit signal - shutting down")
    STOP = True


def decode_bdf(buf):
    """Decode BDF format (0x40)."""
    alert_id = (buf[0] & BDFMASK1)
    alert_id |= (buf[1] & BDFMASK1) << 6
    alert_id |= (buf[2] & BDFMASK3) << 12
    
    alert_data = (buf[2] & BDFMASK2) >> 1
    alert_data |= (buf[3] & BDFMASK1) << 5
    
    return alert_id, alert_data


def decode_eif(buf):
    """Decode EIF format (0xC0)."""
    alert_id = (buf[0] & EIFMASK1)
    alert_id |= (buf[1] & EIFMASK3) << 6
    
    alert_data = (buf[1] & EIFMASK2) >> 7
    alert_data |= (buf[2] & 0xFF) << 1
    alert_data |= (buf[3] & EIFMASK5) << 9
    
    return alert_id, alert_data


def decode_adf(buf):
    """Decode ADF format (0xB0)."""
    id_ones = buf[0] & ADFMASK1
    id_tens = buf[1] & ADFMASK1
    data_ones = buf[2] & ADFMASK1
    data_tens = buf[3] & ADFMASK1
    
    alert_id = id_tens * 10 + id_ones
    alert_data = data_tens * 10 + data_ones
    
    return alert_id, alert_data


def decode_alert_packet(buf):
    """Decode a 4-byte ALERT packet."""
    if len(buf) < 4:
        return None, None, None
    
    first_byte = buf[0]
    
    if (first_byte & 0xC0) == 0xC0:
        alert_id, alert_data = decode_eif(buf)
        return alert_id, alert_data, "EIF"
    elif (first_byte & 0xB0) == 0xB0:
        alert_id, alert_data = decode_adf(buf)
        return alert_id, alert_data, "ADF"
    elif (first_byte & 0xC0) == 0x40:
        alert_id, alert_data = decode_bdf(buf)
        return alert_id, alert_data, "BDF"
    
    return None, None, None


def get_station_info(config_file, station_id):
    """Look up station metadata from configuration file."""
    if not config_file.exists():
        return None
    
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('|')
            if len(parts) >= 12:
                try:
                    if int(parts[0]) == station_id:
                        return {
                            'id': int(parts[0]),
                            'nwsli': parts[1],
                            'name': parts[2],
                            'shef': parts[6],
                            'div': float(parts[7]) if parts[7] else 1,
                            'base': float(parts[8]) if parts[8] else 0,
                            'mult': float(parts[9]) if parts[9] else 1,
                            'adder': float(parts[10]) if parts[10] else 0,
                            'report': int(parts[11]) if parts[11] else 0
                        }
                except (ValueError, IndexError):
                    continue
    return None


def convert_raw_to_value(raw_data, station):
    """Convert raw ALERT counts to engineering units."""
    return (((raw_data / station['div']) + station['base']) * station['mult']) + station['adder']


def create_shef_file(station, obs_value, obs_time, hostname):
    """Create SHEF-formatted file for AWIPS."""
    filename = "SHEFIN-{}.{}".format(station['id'], int(obs_time.timestamp()))
    filepath = ALERTHOME / "data" / filename
    
    shef_line = ".A {} {:02d}{:02d}{:02d} Z DH{:02d}{:02d}{:02d}  /{} {:6.2f}:".format(
        station['nwsli'],
        obs_time.year % 100, obs_time.month, obs_time.day,
        obs_time.hour, obs_time.minute, obs_time.second,
        station['shef'], obs_value
    )
    
    content = "{} ALL\n\n".format(PILID)
    content += ":NATIONAL WEATHER SERVICE {}\n".format(WFOSITE)
    content += ":SUPPLEMENTARY ALERT REPORT\n"
    content += ":SENT FROM {}\n\n".format(hostname)
    content += shef_line + "\n\n"
    
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath, shef_line


def send_to_awips(filepath):
    """Copy file to AWIPS incoming directory."""
    dest = Path("/data/Incoming")
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy(str(filepath), str(dest / filepath.name))


def log_rr5(station_id, shef_line):
    """Append to RR5 data file."""
    RR5DATAFILE.parent.mkdir(parents=True, exist_ok=True)
    with open(RR5DATAFILE, 'a') as f:
        f.write("{}: {}\n".format(station_id, shef_line))


def process_alert_data(alert_id, alert_data, config_file, log_rr5_flag):
    """Process decoded ALERT data - convert and send to AWIPS."""
    now = datetime.now(timezone.utc)
    hostname = socket.gethostname().upper()
    
    station = get_station_info(config_file, alert_id)
    
    if station is None:
        log_message("Station {} not in config".format(alert_id))
        return
    
    obs_value = convert_raw_to_value(alert_data, station)
    
    log_message("ALERTDATA: ID={} raw={} value={:.2f}".format(
        alert_id, alert_data, obs_value))
    
    if station['report'] == 1:
        filepath, shef_line = create_shef_file(station, obs_value, now, hostname)
        send_to_awips(filepath)
        log_message("Sent to AWIPS: {}".format(filepath.name))
        
        if log_rr5_flag:
            log_rr5(alert_id, shef_line)


def main():
    global DEBUG, USE_SYSLOG, STOP, PGMNAME
    
    parser = argparse.ArgumentParser(
        description='Monitor ALERT flood gauge serial port'
    )
    parser.add_argument('-p', '--port', type=int, choices=[1, 2, 3, 4], required=True,
                        help='Port number 1-4')
    parser.add_argument('-b', '--baud', type=int, default=DEFAULT_BAUD,
                        help='Baud rate (default: 9600)')
    parser.add_argument('-s', '--serial', type=str,
                        help='Override serial device path')
    parser.add_argument('-f', '--config', type=Path, default=CONFIG_FILE,
                        help='Station configuration file')
    parser.add_argument('-l', '--logdir', type=Path, default=ALERTLOGDIR,
                        help='Log directory')
    parser.add_argument('-r', '--rr5', action='store_true',
                        help='Log data for RR5')
    parser.add_argument('-D', '--debug', action='store_true',
                        help='Debug mode')
    parser.add_argument('-S', '--syslog', action='store_true',
                        help='Log to syslog')
    args = parser.parse_args()
    
    DEBUG = args.debug
    USE_SYSLOG = args.syslog
    
    # Set program name based on port (matches old naming: n1a, n1c, n1d)
    port_letter = chr(ord('a') + args.port - 1)  # 1->a, 2->b, 3->c, 4->d
    PGMNAME = "alertMon_n1{}".format(port_letter)
    
    # Determine serial device
    serial_device = args.serial if args.serial else PORT_MAP[args.port]
    
    # Setup syslog
    if USE_SYSLOG:
        syslog.openlog(PGMNAME, syslog.LOG_PID, syslog.LOG_LOCAL1)
    
    # Setup signal handlers
    signal.signal(signal.SIGHUP, rotate_log)
    signal.signal(signal.SIGUSR1, exit_handler)
    signal.signal(signal.SIGTERM, exit_handler)
    signal.signal(signal.SIGINT, exit_handler)
    
    # Open log file
    open_log_file(args.logdir, PGMNAME)
    log_message("Starting {} - device={} baud={}".format(PGMNAME, serial_device, args.baud))
    
    # Open serial port
    try:
        ser = serial.Serial(
            port=serial_device,
            baudrate=args.baud,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            timeout=1,
            rtscts=True
        )
        log_message("Opened serial port: {}".format(serial_device))
    except serial.SerialException as e:
        log_message("ERROR: Cannot open {}: {}".format(serial_device, e), syslog.LOG_CRIT)
        return 1
    
    # Main loop
    buffer = bytearray()
    
    while not STOP:
        try:
            data = ser.read(4)
            if not data:
                continue
            
            buffer.extend(data)
            
            while len(buffer) >= 4:
                packet = buffer[:4]
                buffer = buffer[4:]
                
                alert_id, alert_data, fmt = decode_alert_packet(packet)
                
                if alert_id is not None and alert_id > 0:
                    log_message("Received {}: ID={} DATA={}".format(fmt, alert_id, alert_data))
                    process_alert_data(alert_id, alert_data, args.config, args.rr5)
                    
        except serial.SerialException as e:
            log_message("Serial error: {}".format(e), syslog.LOG_ERR)
            break
    
    # Cleanup
    log_message("Shutting down")
    ser.close()
    if logfd:
        logfd.close()
    if USE_SYSLOG:
        syslog.closelog()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
