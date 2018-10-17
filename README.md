# DeepSimulator
The first deep learning based Nanopore simulator which can simulate the process of Nanopore sequencing.

Paper: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty223/4962495

Reference:
```
@article{deepsimulator,
    author = {Li, Yu and Han, Renmin and Bi, Chongwei and Li, Mo and Wang, Sheng and Gao, Xin},
    title = {DeepSimulator: a deep simulator for Nanopore sequencing},
    journal = {Bioinformatics},
    volume = {34},
    number = {17},
    pages = {2899-2908},
    year = {2018},
    doi = {10.1093/bioinformatics/bty223},
    URL = {http://dx.doi.org/10.1093/bioinformatics/bty223},
    eprint = {/oup/backfile/content_public/journal/bioinformatics/34/17/10.1093_bioinformatics_bty223/2/bty223.pdf}
}

```

# Usage
## Prerequisites
1. Anaconda 2 or Minoconda 2 (https://conda.io/miniconda.html)
2. Downloaded Albacore (Suppose you downloaded: ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl. https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl (If the link does not work, please open an issue))

## Build a virtual environment for Albacore basecalling
```
<!-- This the python version should be the same as the Albacore you downloaded -->
conda create -n basecall python=3.6
<!-- Install the Albacore using pip -->
<!-- Make sure the .whl file is in the current directory -->
source activate basecall
pip install ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl
source deactivate
```

## Download the DeepSimulator package
```
git clone https://github.com/lykaust15/DeepSimulator.git
cd ./DeepSimulator/
```

## Build a virtual environment for the context-independent pore model
```
<!-- This step may take up to 8 mins since it would install all the dependencies -->
conda env create -f environment.yml
```

## Run a test
```
./deep_simulator.sh example/human_single.fasta 10
```
<!-- Remember to save the fast5 folder since it is consider to be a temp folder and overwritten every time you run the main.sh file -->

## Control the behavior of DeepSimulator
One can control the behavior of DeepSimulator, including the number of reads generated or the accuracy, *etc.*, by changing the ```main.sh``` file. Detailed descriptions of the parameters in ```main.sh``` file can be refered to Section S4 in [Supplementary material of DeepSimulator](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/34/17/10.1093_bioinformatics_bty223/2/bty223_supplemental_materials.pdf?Expires=2147483647&Signature=v5FSUbbU4eVfQdIo3H3Xrsq6CFVh8azonSxGg1WaAL35maQ0zqIPzdRPTTGUUhlzkLYBnU3Fi4G1DRcXc5YDD4Ea~8ic56zpjBNWQ4qqSZabjH9XwTFyPTbh6IKaHkULi9zKfSl02MxxXfqEJ0Xi72AKRv0Up4j~baWrfyUEKYtEVkzJbpbyAsnZhvPh2WSbFXyPhRhBn~fH9XfrEO9hbMQPSrT9di2Ho85ZBbZ2xS0P~J8sZyi91ulXQHfYSH5rbdaTNAAxRCVbLQUi3iKbJFE5Bl0de66w0mdjfgIZGqrBY9uPoXwW4MYf6H7OwmVXnDc-sjoe73UxO4xHRF17Ag__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

## Advanced
The above example only uses CPU, which would take years to train a model. To accelerate the training process and take advantage the computational power of GPU, users can consider using the GPU version of Tensorflow. User should make sure the following dependencies are installed correctly before running the training code on a workstation with GPU card.

1. CUDA (http://docs.nvidia.com/cuda/cuda-installation-guide-linux/#axzz4VZnqTJ2A)
2. cuDNN (https://developer.nvidia.com/cudnn)
3. Tensorflow-gpu (https://www.tensorflow.org/install/install_linux)

Users can refer to the Tensorflow website (https://www.tensorflow.org/) for more detailed instruction of setting up the environment.

