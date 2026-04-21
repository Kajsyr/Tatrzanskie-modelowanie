include { DATASET } from './modules/dataset.nf'
include { FEATURES_EXTRACTION } from './modules/features.nf'


workflow {
    dataset_ch=DATASET()
    FEATURES_EXTRACTION(dataset_ch.species_ch, dataset_ch.dem_ch)
}