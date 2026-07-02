# config.R
proj_dir  <- file.path("C:/Users/rebec/OneDrive - University College London/Dx test variation/Analysis/Diagnostic waiting times/DxWaitingTimes2")

data_dir  <- file.path(proj_dir, "base_data")
data_out  <- file.path(proj_dir, "data")
output    <- file.path(proj_dir, "output")

# # Shared analysis settings
test_order <- c("Gastroscopy","MRI",
                "Colonoscopy","CT",
                "Flexible Sigmoidoscopy","Non Obstetric Ultrasound",
                "Cystoscopy","Echocardiography")

imaging_tests   <- c("MRI", "CT", "Non Obstetric Ultrasound", "Echocardiography")
endoscopy_tests <- c("Gastroscopy", "Colonoscopy", "Flexible Sigmoidoscopy", "Cystoscopy")

test_order_tables <- c("Gastroscopy", "Colonoscopy", "Flexible Sigmoidoscopy", "Cystoscopy",
                       "MRI", "CT", "Non Obstetric Ultrasound", "Echocardiography")


# # Shared analysis settings - test runs
# test_order <- c("MRI","Gastroscopy")
# imaging_tests   <- c("MRI")
# endoscopy_tests <- c("Gastroscopy")
