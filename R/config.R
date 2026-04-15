# config.R
proj_dir  <- file.path("C:/Users/rebec/OneDrive - University College London/Dx test variation/Analysis/Diagnostic waiting times/DxWaitingTimes2")

data_dir  <- file.path(proj_dir, "base_data")
data_out  <- file.path(proj_dir, "data")
output    <- file.path(proj_dir, "output")

# Shared analysis settings
test_order <- c("MRI","Gastroscopy","CT","Flexible Sigmoidoscopy",
                "Non Obstetric Ultrasound","Colonoscopy","Echocardiography","Cystoscopy")
imaging_tests   <- c("MRI", "CT", "Non Obstetric Ultrasound", "Echocardiography")
endoscopy_tests <- c("Gastroscopy", "Colonoscopy", "Flexible Sigmoidoscopy", "Cystoscopy")