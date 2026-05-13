include { DATASET } from './modules/dataset.nf'
include { BIO_CLIM } from './modules/bioclim.nf'
include { FEATURES } from './modules/features.nf'
include { GEE} from './modules/gee.nf'


workflow {
    dataset_ch=DATASET()
    FEATURES(dataset_ch.full_table,
            dataset_ch.pseudo_absences,
            BIO_CLIM(),
            GEE())
}