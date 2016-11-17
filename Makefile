all:
	( cd src; make )
	( cd cli; make cli )

clean:
	( cd src; make clean )
	( cd cli; make clean )

