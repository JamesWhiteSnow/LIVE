# LIVE

In this project, we open-source the source code and data sets of our Learnable Monotonic Vertex Embedding (LIVE) approach for exact subgraph matching.

## Getting Started

### Dependencies
1. The codes require the following dependences:

* A modern C++ compiler compliant with the C++17 standard (gcc/g++ >= 12.2)
* CMake (>= 3.28)

2. Under the Offline directory, execute the following conda commands to configure the Python environment.

```
conda create --name <new_environment_name> --file requirements.txt
conda activate <new_environment_name>
```

### Offline Process
1. Turn into the Offline directory, execute the following command to train the embedding model.

```
python model_training.py
```

|Parameter|Default Value|Description|
|:----:|:----:|:----:|
|-n|../Dataset|Dataset Path|
|-d|2|Embedding Dimension|
|-e|1000|Epochs|
|-b|4096|Batch Size|
|-l|0.01|Learning Rate|

### Online Process
1. Turn into the Online diectory, execute the following command to build the project.

```
mkdir build
cd build
cmake ..
make
```

2. Return to the root directory, and execute the following command to run a quick start example.

```
./build/LIVE
```

|Parameter|Default Value|Description|
|:----:|:----:|:----:|
|-d|../Dataset/|Dataset Path|
|-q|../Dataset/query_graph.graph|Query Graph Path|
|-a||1000|Alpha Value|
|-b|0.01|Beta Value|
|-e|2|Embedding Dimension|
|-k|2|Hop Number|

