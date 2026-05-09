include { DATASET } from './modules/dataset.nf'
include { BIO_CLIMATIC } from './modules/features.nf'
include { GEE} from './modules/gee.nf'


workflow {
    dataset_ch=DATASET()
    GEE(dataset_ch.full_table)
    //BIO_CLIMATIC(dataset_ch.species_ch, dataset_ch.dem_ch)
}