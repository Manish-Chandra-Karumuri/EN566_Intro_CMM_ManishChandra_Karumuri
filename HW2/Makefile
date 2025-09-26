THETA = 45
WIDTH = 1000

golf:
	python golf.py --plot=$(THETA) 

carbon:
	python carbon.py --plot=$(WIDTH)

.PHONY: all clean

all:
	python golf.py --plot=$(THETA)
	python carbon.py --plot=$(WIDTH)

clean:
	rm -rf *.png