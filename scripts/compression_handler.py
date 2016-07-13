## Compression Handler
##
## Functionality to handle the compression and decompression of images

from os import path
import subprocess
from sys import argv


class DataFile:
    """Data file object with compression-handling functionality"""

    def __init__( self, path_string=None ):
        self.file_path = path_string
        self.extn = None
        self.compress_state = None
        self.output_path = None
	
        # Test for compression by examining the file extension, if any:
        if self.file_path != None:
            (self.output_path, self.extn ) = path.splitext( self.file_path )
        
        # Test for recognized compression states:
        if self.extn == '.fz':
            self.compress_state = 'fzip'
        elif self.extn == '.gz':
            self.compress_state = 'gzip'
        elif self.extn == '.bz2':
            self.compress_state = 'bzip2'
        elif self.extn not in [ '.fits', '.fts' ]:
            self.compress_state = 'unrecognized'
    
    def funzip( self ):
        """Method to uncompress an fzip'd file"""
        
        if path.isfile('/Applications/funpack') == True:
            funpack = '/Applications/funpack'
        else:
            funpack = 'funpack'
        
        if self.compress_state == 'fzip':
            coutput = subprocess.check_call( [ funpack + ' ' + str(self.file_path) ], \
                    stderr=subprocess.STDOUT, shell=True )
        else:
            coutput = 'INFO: File not fzipped, no decompression performed'

def uncompress( file_path, out_file_path=None ):
    """Function to uncompress a single file"""
    data_file = DataFile( path_string=file_path )
    data_file.funzip()

if __name__ == '__main__':
    
    if len(argv) > 1:
      	file_path = argv[1]
    else:
      	file_path = raw_input('Please enter file_path: ')
    uncompress( file_path )
    
    