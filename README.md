# CAncer-bioMarker-Prediction-Pipeline - CAMPP2  #
<br/>

### Installation instructions

#### Create and activate conda environment with R and devtools in your project directory
`conda create --prefix env -c conda-forge r-base=4.1.3 r-devtools`
`conda activate ./env`


#### Install CAMPP2 from Github
open R console: <br/>
`R` <br/>


install dependency for devtools <br/>
`install.packages("usethis")`

load devtools library <br/>
`library(devtools)`

install CAMPP2 from private Github repository (using your personal token and commit ID)
`devtools::install_github(repo = "ELELAB/CAMPP2",ref="668e7b6b8a7b69a711a961b0c96d75dd3afb7d30",auth_token="<your personal token to github>")`

### Example run
Default settings can be executed in R using 
```
library(CAMPP2)
library(edgeR)

runCampp2(batches=c("tumor_stage","tumor_stage"),prefix="test_CAMPP2", data1=dataset1, data2=dataset2, metadata1=metadata1,metadata2=metadata2, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"))

```
<br/>

For more details, see help page of the main function `runCampp2`.

