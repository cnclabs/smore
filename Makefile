all:
	( cd src; make )
	( cd cli; make all )

clean:
	( cd src; make clean )
	( cd cli; make clean )

