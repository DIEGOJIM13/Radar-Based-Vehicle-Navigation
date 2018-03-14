import serial
import subprocess
from time import gmtime, strftime

# Usage 
# Run: python -m serial.tools.list_ports
# 
# $GPRMC,hhmmss.ss,A,llll.ll,a,yyyyy.yy,a,x.x,x.x,ddmmyy,x.x,a*hh
#
# 1    hhmmss.ss 	= UTC of position fix
# 2    A 			= Data status (V=navigation receiver warning)
# 3    llll.ll 		= Latitude of fix
# 4    a 			= N or S
# 5    yyyyy.yy 	= Longitude of fix
# 6    a 			= E or W
# 7    x.x 			= Speed over ground in knots
# 8    x.x 			= Track made good in degrees True
# 9    ddmmyy 		= UT date
# 10   x.x 			= Magnetic variation degrees (Easterly var. subtracts from true course)
# 11   a 			= E or W
# 12   hh 			= Checksum

def main():
	output = subprocess.check_output(['python','-m', 'serial.tools.list_ports'])
	outList = output.split('\n')
	outList = outList[0:-1]
	print(outList)
	for serialName in outList:
		print(serialName)  
		ser = serial.Serial(timeout=2)
		ser.baudrate = 115200
		ser.port = serialName
		try:
			ser.open()
		except Exception:
			print('Exception in opening')
		line = ser.readline()        
		if (len(line) is not 0):
			fileName = strftime("%m-%d-%H-%M-%S", gmtime())  + '.txt'
			if (ser.is_open):
				while (True):
					input = ser.readline()
					input.replace('\n', '')                    
					if ('Start' in input):
    						print('Start')
    						f = open(fileName,'a')
    						f.write('Start\n')
    						f.close()   
					elif ('End' in input):
    						print('End')
    						f = open(fileName,'a')
    						f.write('End\n')
    						f.close()                              
					elif ('$GPRMC' in input):
						arr = [x.strip() for x in input.split(',')]
						if (len(arr) == 13):
        						line = arr[0] + ',' + arr[1] + ',' + arr[2] + ',' + arr[3] + ',' + arr[4] + ',' + arr[5] + ',' + arr[6] + ',' + arr[7] + ',' + arr[8] + ',' + arr[9] + ',' + arr[10] + ',' + arr[11] + ',' + arr[12]
        						print(line)
        						f = open(fileName,'a')
        						f.write(line.replace('\n', '') + '\n')
        						f.close()                        
			else:
				print('Could not open that serial port')
		ser.close()
			
	print('Done')
	
if __name__ == '__main__':
	main()
