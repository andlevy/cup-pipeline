from datetime import datetime, timedelta
import subprocess
import os
import argparse


# This code is activated via a crontab
# which can be seen using sudo cat /var/spool/cron/crontabs/skw24 on pearlai (139.184.160.16)
# or sudo cat /var/spool/cron/crontabs/skw24aza on linuxbox (139.184.20.69)
# The job is run every night
# crontab -e

def is_directory_expired(dir_age: datetime):
    expiration_date = datetime.now() - timedelta(days=6)
    return dir_age < expiration_date


def find_age_of_directory(patient_dir):
    stat_info = os.stat(patient_dir)
    modification_timestamp = stat_info.st_mtime
    modification_datetime = datetime.fromtimestamp(modification_timestamp)
    return modification_datetime


def find_expired_paths(patient_directories):
    expired_paths = []
    for patient_dir in patient_directories:
        dir_age = find_age_of_directory(patient_dir)
        if is_directory_expired(dir_age):
            expired_paths.append(patient_dir)
    return expired_paths


def remove_expired_directories(data_dir):

    """Identify those directories in the data directory which have not been marked as keep
    and delete them if they are more than 6 days old"""
    if data_dir[-1] == "/":
        data_dir = data_dir[0:-1]

    def data(file):
        return f"{data_dir}/{file}"

    patient_directories = [data(p) for p in os.listdir(data_dir) if os.path.isdir(data(p))]
    expired_dirs = find_expired_paths(patient_directories)
    for directory in expired_dirs:
        subprocess.run(['rm', '-r', directory])
    expired_dirs_string = "\n".join(expired_dirs) + '\n'
    message = ("This cron job was last run on {date}. "
               "It inspected {data_dir} for old files. "
               "The following directories were removed: {lost_dirs}".format(date=datetime.now().strftime("%d_%m_%Y"),
                                                                            data_dir=data_dir,
                                                                            lost_dirs=expired_dirs_string))
    keep_txt = data('keep.txt')
    if not os.path.isfile(keep_txt):
        with open(keep_txt, 'w') as f:
            f.write(message)
    else:
        with open(keep_txt, 'a') as f:
            f.write(message)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Remove patient data that is more than 7 days old")
    parser.add_argument('--data_dir', type=str, help='directory that stores patient data')

    args = parser.parse_args()
    remove_expired_directories(args.data_dir)
