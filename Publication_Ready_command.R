library(here)
library(rmarkdown)

Model<-c("Model_1","Model_2"); i=1

for(i in 1:length(Model))
{
  # Publication Ready ----
  render(input = here("Publication_Ready.Rmd"),
         output_format = "github_document",
         output_file = "Publication_Ready",
         output_dir = here("Publication_Ready",Model[i]),
         params=list("Model_Path"=Model[i]))

}
