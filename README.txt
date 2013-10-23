####################################################
HOWTO
####################################################



# Things you need:

1. Enthought Python Distribution (or scipy)

2. Copy python/TerminalController.py into your python sys.path.
Do this by:
	$ python
	>> import sys, pprint
	>> pprint.pprint( sys.path )

	and pick a directory you like from that list

	$ sudo cp python/TerminalController.py /directory/you/liked/


3. The C++ BOOST library from http://boost.org -- particularly boost/foreach.hpp
