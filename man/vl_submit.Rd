% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vl_submit.R
\name{vl_submit}
\alias{vl_submit}
\title{Submit Commands to a Server with Resource Management}
\usage{
vl_submit(
  cmd,
  cores,
  mem = 16,
  time = "08:00:00",
  job.name,
  logs = "db/logs",
  overwrite = FALSE,
  execute = TRUE
)
}
\arguments{
\item{cmd}{A data.table containing the commands to execute, with columns:
- `file.type`: Type of output file.
- `path`: Path to the output file.
- `cmd`: Shell command to execute.
- `cores`: (Optional) The number of CPU cores to allocate to the job, specified by the first value.
- `job.name`: (Optional) The name of the job, specified by the first value in this column.}

\item{cores}{Number of CPU cores to allocate. If not specified, the first value of the optional
'core' column will be used. Otherwise, default= 8.}

\item{mem}{Memory allocation in GB. Default= 16.}

\item{time}{Maximum runtime for the job (HH:MM:SS). Default= '08:00:00'.}

\item{job.name}{Name of the job. If not specified, the first value of the optional 'job.name'
column will be used. otherwise, default= "myJob".}

\item{logs}{Directory for log files. Default: 'db/logs'.}

\item{overwrite}{If set to TRUE, overwrites existing output files. Default= FALSE.}

\item{execute}{If set to FALSE, the command is returned and is not submitted to the cluster.
Default= TRUE.}
}
\description{
Submits shell commands to a server using resource management (e.g., LSF). Handles job submission,
resource allocation, and optional directory creation for output files.
}
\examples{
# Submit a command to the server
cmd <- data.table(
  file.type = "bam",
  path = "output.bam",
  cmd = "echo 'Hello world'",
  job.name= "test"
)
vl_submit(cmd)

# Return the constructed command without submitting
vl_submit(cmd, execute = FALSE)

}
