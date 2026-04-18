include { DATASET } from './modules/dataset.nf'
//include { FEATURES } from './modules/features.nf'


workflow {
    DATASET()
    //FEATURES()
}