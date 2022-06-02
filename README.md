# CAncer-bioMarker-Prediction-Pipeline - CAMPP2  #
<br/>

### Installation instructions

#### Create conda environment with R and devtools
`conda create --prefix env -c conda-forge r-base=4.1.3 r-devtools`

#### Install CAMPP2 from Github
open R console: <br/>
`R` <br/>


install dependency for devtools <br/>
`install.packages("usethis")`

load devtools library <br/>
`library(devtools)`

install CAMPP2 from private Github repository (using your personal token and commit ID)
`devtools::install_github(repo = "ELELAB/CAMPP2",ref="81bd7b3207ea49a1d5331e7d28e2201f66d96212",auth_token="<your personal token to github>")`

### Example run
Default settings can be executed in R using `campp2_exampleRun.R`. 
<br/>

For more details, see help page of the main function `runCampp2`.

