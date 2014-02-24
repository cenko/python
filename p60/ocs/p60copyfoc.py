#! /usr/bin/env python

import sys, socket, time, os, re, string, shutil, ephem
import telnetlib, getopt, oscar_session, mx.DateTime
import smtplib
from email.MIMEText import MIMEText
from types import *
from math import *
from iqutils import *

rawfiles = {}
dir = os.listdir('.')
for file in dir:
    rawfiles[file] = 1

while (time.localtime()[3] > 15 or time.localtime()[3] < 8) :

    dir = os.listdir('.')
    for file in dir:
        if re.search('20\d*r.fits', file) and (not rawfiles.has_key(file)):
            time.sleep(30)
            [imgtype] = get_head(file, ['IMGTYPE'])
            if (imgtype == 'SAOFOCUS'):
                os.system('cp %s ../hipfocus' % file)
            rawfiles[file] = 1
    time.sleep(60)

