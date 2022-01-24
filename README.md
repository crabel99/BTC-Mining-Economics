# Bitcoin Marginal Utility

This repository presents a simple model for evaluating the marginal utility of bitcoin.
The analysis work is done in RStudio.
The onchain bitcoin data is parsed from the bitcoin blockchain using <https://github.com/blkchain/blkchain>.
The parsed bitcoin blockchain data is stored in a PostgreSQL database.

To be able to use this code you will need to install:

- [PostgreSQL](https://www.postgresql.org/download/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/) or [R](https://www.r-project.org/)
- A blockchain parser e.g. [Blockchain Postgres Import](https://github.com/blkchain/blkchain)

If you are using RStudio , you will need to enter `renv::restore()` in the console when you first open the project.
Once you have the environment set up and the database populated, you will need to open `Bitcoin Marginal Utility.R` and execute it from source.
This will run the entire file and provide you wil the graphical outputs seen in the `*.pdf`'s in this repository.

The Bitcoin price history was downloaded from the [St. Louis Federal Reserve Bank](https://fred.stlouisfed.org/series/CBBTCUSD) and only goes back to 2014.
If you know of a more comprehensive data source. please open an issue so that I can update the model.

This project is covered under [Creative Commons 4.0-BY](https://creativecommons.org/licenses/by/4.0/legalcode).
Please feel free to do with this as you will, but please remember to cite this work appropriately.
