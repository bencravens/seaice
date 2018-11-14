"""script file to plot grabbed variables over time"""
import grab
import mystats

def main():
	grab.get("/media/windowsshare","aice","u-au866","20060401-20060501")

if __name__=='__main__':
	main()
