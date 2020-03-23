'''This code is used only to retrieve the urls of the zip files of the faers, it also install the necesssary packages if missed'''

import sys
import subprocess

try:
	from bs4 import BeautifulSoup
	import re
	import requests
except ImportError:
	subprocess.call([sys.executable, "-m", "pip", "install", 'bs4'])
	subprocess.call([sys.executable, "-m", "pip", "install", 're'])
	subprocess.call([sys.executable, "-m", "pip", "install", 'request'])
finally:
	from bs4 import BeautifulSoup
	import re
	import requests

response = requests.get('https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html')
soup = BeautifulSoup(response.content, 'html.parser')
for link in soup.findAll('a'):
    print(link.get('href'))
