# MutAnt

Current code contains code to convert dbNSFP database into suitable for MutAnt format

## Usage

Download the latest dbNSFP release from the official [site](https://sites.google.com/site/jpopgen/dbNSFP). We recommend to use softgenetics ftp link.

```bash
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.6a.zip
```
Unzip the downloaded folder with any suitable program and locate into folder. For example:
```bash
unzip dbNSFP4.6a.zip && cd dbNSFP4.6a
```

Copy file suitable for chosen genome version, hg38 or hg19 into folder with dbNSFP database.
From the folder run:
```bash
cp ~/MutAnt/generate_aggregated_database_hg38.sh .
```
Edit dbNSFP version in the downloaded file, default version is set to 4.6a. For, example, open file with `nano` and change version in line 3: `version="4.6a"`

Launch the conversion script
```bash
sh generate_aggregated_database_hg38.sh
```