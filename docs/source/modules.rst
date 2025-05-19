Modules de wisp_light
===========================

.. toctree::
   :maxdepth: 2
   :caption: Modules



wisp_light.build_dataset.refseq
-------------------------------

.. literalinclude:: ../../wisp_light/build_dataset/refseq/download_refseq_from_csv.py
   :language: python
   :linenos:
   :caption: Script download_refseq_from_csv.py

.. literalinclude:: ../../wisp_light/build_dataset/refseq/get_all_taxo_from_NCBI.py
   :language: python
   :linenos:
   :caption: Script get_all_taxo_from_NCBI.py


wisp_light.dataset
--------------------------------

.. automodule:: wisp_light.dataset.refSeqDataset
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: wisp_light.dataset.bactero_set
   :members:
   :undoc-members:
   :show-inheritance:

wisp_light.training
-------------------

Ce sous-package contient les modules liés à la création de bases de données, de modèles et de prédictions.

wisp_light.training.create_database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: wisp_light.training.create_database
   :members:
   :undoc-members:
   :show-inheritance:

wisp_light.training.create_model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: wisp_light.training.create_model
   :members:
   :undoc-members:
   :show-inheritance:

wisp_light.training.create_prediction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: wisp_light.training.create_prediction
   :members:
   :undoc-members:
   :show-inheritance:

wisp_light.training.training_functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: wisp_light.training.training_functions
   :members:
   :undoc-members:
   :show-inheritance:

wisp_light.training.metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: wisp_light.training.metrics
   :members:
   :undoc-members:
   :show-inheritance:


wisp_light.prediction
---------------------

.. automodule:: wisp_light.prediction.predict
   :members:
   :undoc-members:
   :show-inheritance:


Script d'entraînement
~~~~~~~~~~~~~~~~~~~~~

Ce script lance l’entraînement du modèle bactérien, avec suivi via MLflow, gestion des logs, création ou rechargement de la base de données.

.. literalinclude:: ../../wisp_light/training/train_val.py
   :language: python
   :linenos:
   :caption: Script train_val.py

wisp_light.visu
---------------
.. automodule:: wisp_light.visu.plots_tools
   :members:
   :undoc-members:
   :show-inheritance:
