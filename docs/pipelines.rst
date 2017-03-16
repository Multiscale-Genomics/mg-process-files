Pipelines
=========

BED File Indexing
-----------------
.. automodule:: process_bed
   
   This pipeline can process bed files into bigbed and HDF5 index files for web
   use.
   
   Running from the command line
   =============================
   
   Parameters
   ----------
   assembly : str
      Genome assembly ID (e.g. GCA_000001405.22)
   chrom : int
      Location of chrom.size file
   bed_file : str
      Location of input bed file
   h5_file : str
      Location of HDF5 output file
   
   Returns
   -------
   BigBed : file
      BigBed file
   HDF5 : file
      HDF5 index file
   
   Example
   -------
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
   
   ``chrom.size`` file:
   .. code-block:: none
      :linenos:

      1  123000000
      2  50000000
      3  25000000
      4  10000000
      5  5000000
      X  75000000
      Y  12000000

   .. code-block:: none
      :linenos:
      
      runcompss --lang=python /home/compss/mg-process-files/process_bed.py --assembly GCA_000001405.22 --chrom chrom.size --bed_file <data_dir>/expt.bed --h5_file <data_dir>/expt.hdf5

   Process Methods
   ---------------
   .. autoclass:: process_bed.process_bed
      :members:

WIG File Indexing
-----------------
.. automodule:: process_wig
   
   This pipeline can process WIG files into bigbed and HDF5 index files for web
   use.
   
   Running from the command line
   =============================
   
   Parameters
   ----------
   assembly : str
      Genome assembly ID (e.g. GCA_000001405.22)
   chrom : int
      Location of chrom.size file
   wig_file : str
      Location of input wig file
   h5_file : str
      Location of HDF5 output file
   
   Returns
   -------
   BigWig : file
      BigWig File
   HDF5 : file
      HDF5 index file
   
   Example
   -------
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
   
   ``chrom.size`` file:
   .. code-block:: none
      :linenos:

      1  123000000
      2  50000000
      3  25000000
      4  10000000
      5  5000000
      X  75000000
      Y  12000000

   .. code-block:: none
      :linenos:
      
      runcompss --lang=python /home/compss/mg-process-files/process_wig.py --assembly GCA_000001405.22 --chrom chrom.size --wig_file <data_dir>/expt.wig --h5_file <data_dir>/expt.hdf5
   
   Process Methods
   ---------------
   .. autoclass:: process_wig.process_wig
      :members:

3D JSON Indexing
----------------
.. automodule:: process_json_3d
   
   This pipeline processes the £D JSON models that have been generated via
   TADbit into a single HDF5 file that can be used as part of a RESTful API for
   efficient querying and retrieval of the models.
   
   Running from the command line
   =============================
   
   Parameters
   ----------
   gz_file : str
      Location of the input tar.gz file containing all of the output models and
      data from the TADbit modelling stage.
   
   Returns
   -------
   HDF5 : file
      HDF5 index file
   
   Example
   -------
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):
   
   .. code-block:: none
      :linenos:
      
      runcompss --lang=python /home/compss/mg-process-files/process_json_3d.py --gz_file <data_dir>/expt.tar.gz

   Process Methods
   ---------------
   .. autoclass:: process_json_3d.process_json_3d
      :members:

