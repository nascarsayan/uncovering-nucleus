# Uncovering Nucleus

This project is the python implementation of the paper "Uncovering the Nucleus of Social Networks by Braulio Dumba and Zhi-li Zhang"

## Installing dependencies

* Requires python3 and pip3
* Install the required python packages for the project using
```bash
pip3 install -r requirements.txt
```

## Usage

```bash
python3 init.py ['/path/to/edgelist/file']
```

For example,

```bash
python3 init.py './dataset/ca-AstroPh.txt'
```

This will output the Nuclear-Index graph in the folder './out/dataset/ca-AstroPh.txt/ni'