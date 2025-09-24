.. DiffusionEquationCode documentation master file, created by
   sphinx-quickstart on Fri Apr 14 14:14:54 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Diffusion Equation's documentation
=================================================

The software is an open source new generation of room acoustics software for research, education and industry in acoustics.

The software is based on the Acoustics Diffusion Equation Method for modeling of sound behaviour in complex geometrical spaces.
The Diffusion Equation Method is used to understand the acoustics properties of the room and to obtain spatial distribution of acoustics energy over time in specific rooms.
The software is developed with two numerical methods: the Finite Different Method and the Finite Volume Method.
Each method is distributed with its own Python code. The main application of the method is room and building acoustics, but in the future it could also lead to more applications.


.. toctree::
   :maxdepth: 2
   :caption: Software Presentation:

   Overview.md

.. toctree::
   :maxdepth: 2
   :caption: Software Use:

   Finite Different Method Use.md
   Finite Volume Method Use.md
   Auralization Use.md

.. toctree::
   :maxdepth: 2
   :caption: Software Theory:

   DocumentationFDM.md
   DocumentationFVM.md
   DocumentationAuralization.md

.. toctree::
   :maxdepth: 6
   :caption: API reference:

   acousticDE/FDMfunctions
   acousticDE/FVMfunctions
   acousticDE/ReverberationFunctions
   acousticDE/Auralizationfunctions
   
.. toctree::
   :maxdepth: 2
   :caption: Tutorials:

   Finite Different Method Demo Tutorial.md
   Finite Volume Method Demo Tutorial.md

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

