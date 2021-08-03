# soma-germline-timing
R script to analyze C. elegans somatic and germline timings as described in the paper 'Neuronal perception of the social environment generates an inherited memory that controls the development and generation time of C. elegans' 

CellProfiler 3 project files with pipelines for measuring fluorescence in C elegans embryos ('Embryofluorescence.cpproj') or for measuring the length of L1 larvae. 

To measure length, it is first necessary to produce a training set for use by the UntangleWorms module of the CellProfiler WormToolbox package. The pipelines 'L1length_makebinaryfortraining.cpproj' and 'L1length_trainingpipeline.cpproj' are used to create this training set. The training set can be produced using the images for measurement. Once the training set is created the worms can be measured with 'L1length_measureworms.cpproj'. These pipelines are easy to adapt to worms of any developmental stage.
