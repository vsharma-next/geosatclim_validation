# Validation framework for geosatclim

## Purpose

a simple python application to perform validation of satellite data processing using point measurements from different station networks. 

By simple, I mean that I will first start by procedural style scripts to manually run a data pipeline: 

1. Obtain the raw data hopefully programatically or documented steps 
2. Convert raw station data to 'analysis ready data'. We will get into analysis ready data and its meaning later. 
3. calculate simple statistical scores. 

## Station Networks

The stations networks are

- the BSRN network,
- the GEBA network and,
- the SwissMetNet.

These come in multiple forms:

- netcdfs, csv files, etc. 
- different temporal granularity
- gaps

### BSRN

- Metadata = bsrn_stations.tab

```bash
wget -O bsrn_stations.tab https://dataportals.pangaea.de/bsrn/stations?format=textfile
```

lftp -e "set ssl:verify-certificate no; mirror --parallel=8 / .; quit" -u username,password ftp://ftp.bsrn.awi.de

### GEBA

- Metadata = GEBA_2024-12-03_13-53-10_metadata.csv

### SwissMetNet

/TODO

## Satellite Processing

The satellite processing comes in the form of NETCDF files where variables ( surface incoming radiation, albedo etc ) 3D arrays ( x,y,t ) are present.

Important to note is that one file represents one time window and there are multiple such files. Thus, due to the file organization, data is naturally 'sliced' along the time.

## Python Environment

We will ONLY use conda as package / environment manager.

For netcdf handeling - xarray
For station data handling - pandas
For plotting - matplotlib, seaborn, plotly ( no one framework chosen yet )
