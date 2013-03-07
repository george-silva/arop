
MODLIST = source

all:
	@for i in $(MODLIST); do \
	(cd $$i; make) \
	done

install:
	@for i in $(MODLIST); do \
	(cd $$i; make install) \
	done

clean:
	@for i in $(MODLIST); do \
	(cd $$i; make clean > /dev/null) \
	done
