# cfDNAmethylome
Workflow for Twist targeted methylation sequencing of cfDNA 

md5 (Message Digest Algorithm 5)
==================================
+ md5 creates md5 hash
+ To get md5 hash for multiple files in a directory with .fa files, use a wildcard to match all files as follows:

		md5 /check_this_dir/* > files.md5 

md5sum
=============
+ use the md5sum command to verify MD5 hashes for multiple files

		$ md5sum * | tee md5sum.txt 

md5sum compute the hash values of the files in the current directory `*` and
`tee`is used to simultaneously display the outuput on stdout and write the otutput to `md5sum.txt` file


+ run md5sum -c command with the text file as argument

		$ md5sum -c mdsum.txt
		
-`c` option of `md5sum` can accept a text file as an input
