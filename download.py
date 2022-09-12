# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 23:30:28 2020

@author: User
"""

import sys
import urllib
import urllib.request
import fnmatch
import lxml.html
import webbrowser
import time

def url_lister(url):
    urls = []
    connection = urllib.request.urlopen(url)
    dom =  lxml.html.fromstring(connection.read())
    for link in dom.xpath('//a/@href'):
        urls.append(link)
    return urls


def progress_hook(out):
    """Return a progress hook function, suitable for passing to
    urllib.retrieve, that writes to the file object *out*."""
    def it(n, bs, ts):
        got = n * bs
        if ts < 0:
            outof = ''
        else:
            # On the last block n*bs can exceed ts, so we clamp it
            # to avoid awkward questions.
            got = min(got, ts)
            outof = '/%d [%d%%]' % (ts, 100 * got // ts)
        out.write("\r  %d%s" % (got, outof))
        out.flush()
    return it
folder_url = "https://e4ftl01.cr.usgs.gov/MOTA/MCD19A2.006/"
folder_urls = url_lister(folder_url)
folderid="2013.*.*/"
folder_num = [filename for filename in fnmatch.filter(folder_urls, folderid)]
folder_num
for j in range (0,len(folder_num)):

    url = "https://e4ftl01.cr.usgs.gov/MOTA/MCD19A2.006/"+folder_num[j]
    urls = url_lister(url)

    filetype = "*.hdf"
   # fileid="*.h26v06.*.*.hdf"   #BD
    fileid="*.h26v06.*.*.hdf"   #BD
   # fileid="*.h11v04.*.*.hdf"    #Pittsburgh
    file_list = [filename for filename in fnmatch.filter(urls, fileid)]
    print(file_list)
    try:
        for i in range(0,len(file_list)):
        #for i in range(0,2):

            New_url = url+file_list[i]
            hdf_file = url.split('/')[-1]
            #hdf_file = file_list[i]
            #hdf_file = "MCD19A2.A2000057.h35v09.006.2018013034447.hdf"
            sys.stdout.write(hdf_file + '\n')
            #urllib.request.urlretrieve(url, filename=hdf_file, allow_redirects=True
            #                       reporthook=progress_hook(sys.stdout))
            print(New_url)
            webbrowser.open(New_url, autoraise=True)
            sys.stdout.write('\n')
            sys.stdout.flush()
            time.sleep(5)
    except:
        print("Couldnot download file:")
        print(file_list)
