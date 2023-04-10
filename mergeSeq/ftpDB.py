"""
能运行，但是非常的慢
"""
import os
from datetime import datetime
import ftplib
import time


date = datetime.now().strftime("%Y%m%d")

# create a folder named "bacteria.16SrRNA-{date}"

folder_name = f"bacteria.16SrRNA-{date}"
if not os.path.exists(folder_name):
    os.mkdir(folder_name)

# download the file from the given url and save it in the folder
# connect to the FTP server
ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
ftp.login()
print("Successfully connected to FTP server")
ftp.cwd("/refseq/TargetedLoci/Bacteria")
ftp.maxline = 1048576 # set buffer size to 1 MB
filename = "bacteria.16SrRNA.fna.gz"
local_filename = os.path.join(folder_name, filename)

ftp.retrbinary(f"RETR {filename}", open(local_filename, "wb").write)
ftp.quit()
