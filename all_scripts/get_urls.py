'''This code is used only to retrieve the urls of the zip files of the faers, it also install the necesssary packages if missed'''

from bs4 import BeautifulSoup
import requests

response = requests.get('https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html')
soup = BeautifulSoup(response.content, 'html.parser')
for link in soup.findAll('a'):
    print(link.get('href'))
