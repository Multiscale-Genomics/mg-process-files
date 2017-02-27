Tools to index genomic files
============================================

.. automodule:: tool
   
   Bed Indexer
   ----------------
   .. autoclass:: tool.bed_indexer.bedIndexerTool
      :members:
      
      .. automethod:: tool.bed_indexer.bedIndexerTool.bedsort(self, file_bed, file_sorted_bed)
      .. automethod:: tool.bed_indexer.bedIndexerTool.bed2bigbed(self, file_bed, file_chrom, file_bb)
      .. automethod:: tool.bed_indexer.bedIndexerTool.bed2hdf5(self, file_id, assembly, file_sorted_bed, file_hdf5)
