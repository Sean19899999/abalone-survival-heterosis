# Read the revision workbook or the equivalent public CSV snapshot.
read_sheet <- function(path, sheet, skip=0, ...) {
  if (file.exists(path)) return(readxl::read_excel(path, sheet=sheet, skip=skip, ...))
  csv <- file.path(dirname(path), 'data', paste0(sheet, '.csv'))
  if (!file.exists(csv)) stop(paste('Missing exported sheet:', basename(csv)))
  readr::read_csv(csv, skip=skip, show_col_types=FALSE, progress=FALSE, ...)
}
