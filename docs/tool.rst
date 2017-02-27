Tools to index genomic files
============================================

.. automodule:: tool
   
   Bed Indexer
   ----------------
   .. autoclass:: tool.bed_indexer.bedIndexerTool
      :members:
      
      .. automethod:: bedsort(self, file_bed, file_sorted_bed)
      .. automethod:: bed2bigbed(self, file_bed, file_chrom, file_bb)
      .. automethod:: bed2hdf5(self, file_id, assembly, file_sorted_bed, file_hdf5)
