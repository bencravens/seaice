"""main script to use icedata module to analyze data"""
import argparse

def main(work_dir, optional_int):
    print "your work dir is {}".format(work_dir)
    print "your integer is {}".format(optional_int)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='grab data, process it, and make plots')
    parser.add_argument('--work_dir', help="working directory to find files", required=True)
    parser.add_argument('--optional_int', help="optional integer", type=int, default=42)
    args = parser.parse_args()
    print "your args are {} {}".format(args.work_dir,args.optional_int)
    main(args.work_dir,args.optional_int)
