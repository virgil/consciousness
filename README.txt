This is the codebase for computing effective information (ei), the minimum information partition (MIP), integrated information (phi), and other measures related to the _Information Integration Theory_ of Consciousness (IITC).

--------

<img src="http://dl.dropbox.com/u/3308162/glassware.jpg" align="right" alt="" border="0" title="glassware image" width="219" height="234">

==Popular Reading==
  # [https://www.youtube.com/watch?v=OdtT0S0GNjE Howto for developers]
  # *[http://www.nytimes.com/2010/09/21/science/21consciousness.html?_r=2&pagewanted=all New York Times - Sizing Up Consciousness by Its Bits]*
  # *[http://www.scientificamerican.com/article.cfm?id=a-theory-of-consciousness Scientific American Mind: A "Complex" Theory of Consciousness]*
  # [http://www.spectrum.ieee.org/sing_koch Video: Teaching Machines to Watch Blade Runner]
  # [http://www.spectrum.ieee.org/jun08/6278 Can Machines Be Conscious?]
  # [http://www.spectrum.ieee.org/jun08/6315 Consciousness as Integrated Information]

==Directly Relevant Academic Papers==
  # *[http://www.biolbull.org/cgi/content/abstract/215/3/216 Consciousness as Integrated Information: a Provisional Manifesto]*
  # [http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000091 Integrated Information in Discrete Dynamical Systems: Motivation and Theoretical Framework]
  # [http://www.theassc.org/documents/practical_measures_of_integrated_information_for_stationary_systems Practical measures of integrated information for stationary systems]
  # *[http://arxiv.org/abs/1401.0978 A Principled Infotheoretic φ-like Measure]*


==Relationship between IIT and other fields==
 # [http://ntp.neuroscience.wisc.edu/faculty/fac-art/tononicon&anesth.pdf Neuroscience and IIT]

==Related Papers==
  # [http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=331407 Measuring information integration]
  # [http://www.biomedcentral.com/1471-2202/5/42 An information integration theory of consciousness]
  # [http://www.sciencemag.org/cgi/content/abstract/sci;282/5395/1846 Consciousness and Complexity]
  # [http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000462 Qualia: The Geometry of Integrated Information]

------

# About
===============
This is code underlying my [http://thesis.library.caltech.edu/8041/](PhD thesis).
Aside from complexity work related to the [http://arxiv.org/abs/1401.0978](psi measure),
it is *by far* the fastest implementation of [http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000091](Balduzzi-Tononi-2008 phi).  The major restrictions about the code are:

* It only works for DETERMINISTIC networks.
* It only works for networks with <= 32 nodes.



# Requirements
===============
1. The C++ BOOST library from http://boost.org -- particularly boost/foreach.hpp

If on a mac, you can install these using 'macports' or 'homebrew'.


# Compiling
===============
Run the command 'recompile', which generates an executable 'consciousness'.  I.e.,

$ ./recompile
$ ./consciousness filename.txt

Computes the measures for the system specified in filename.txt





# Some notes
=======================

The directory "e/" contains example systems to compute the phi/psi.  For example, to compute the measures for the system "transitions/3RN.txt", you'd do:

$ ./consciousness e/transitions/3RN.txt

-----

The directory 'balduzzi_python' contains the original python code from David Balduzzi to compute the phi in the 2008 paper, "Integrated Information in Discrete Dynamical Systems"

-----

The directory 'tests' is a series of simple programs that spit out system diagnostic information.  It's unlikely you'll ever need them.  You can safely ignore this directory.

-----

The directory 'pics' contains two pretty pictures for the Google Code site.  You can ignore this directory.

-----

