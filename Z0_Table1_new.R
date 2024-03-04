library(tidyverse)
library(gt)
library(webshot2)
library(gtExtras)
tabledata<-read.table("C:/Users/rochf/Desktop/veggiemeat/veggiemeat/01-metadata/table1.txt",sep="\t",header=T, nrow=32)
tabledata<-tabledata[order(tabledata$protein.source,tabledata$texture),]
tabledatared<-tabledata[,c(2,1,5,6,7,8)]
gt::gt(tabledatared, "latex"
 )%>%
  tab_footnote(
    footnote = "frozen product",
    locations = cells_body(columns = ID, rows = 25:26)
  )%>%
  tab_footnote(
    footnote = "Protein basis of the examined product. Only pea or soybean protein products were selected for the study.",
    locations = cells_stubhead()
  )%>%
  tab_footnote(
    footnote = "product designation. Products with a minced 'meat' basis (i.e. minced meat, burger, cevapcici, sausages) were additionally classified as 'minced', products immitating pieces of meat or a meat structure (i.e. fillets, steaks, chunks, kebab) were classified as 'fibrous'.",
    locations = cells_row_groups(groups = everything())
  )%>%
  tab_footnote(
    footnote = "days to expiration date (ed) or best before date (bbd) at sampling. In brackets: consume within x days after opening.",
    locations = cells_column_labels(columns = shelf.life)
  )%>%
  tab_footnote(
    footnote = "Recommended cooking time. If label said (e.g.) 2 minutes per side, the recommended cooking time were doubled to 4 minutes for this table. ",
    locations = cells_column_labels(columns = cooking.time)
  )%>%
  
  tab_row_group(
    label = "pea protein with fibrous texture",
    rows = 1:6,
    id="ppft"
  ) %>%
  tab_row_group(
    label = "pea protein with minced texture",
    rows = 7:16,
    id="ppmt"
  ) %>%
  tab_row_group(
    label = "soybean protein with fibrous texture",
    rows = 17:25,
    id="spft"
  ) %>%
  tab_row_group(
    label = "soybean protein with minced texture",
    rows = 25:32,
    id="spmt"
  ) %>%
  cols_label(
    shelf.life="shelf life",
    cooking.time = "cooking time", 
    no..Of.ingredients = "no. of ingredients",
    additional.labelling = "additional labelling"
  )%>%
  tab_options(
    table.font.name = "Optima",
    table.font.color = "black",
    table.border.top.style = "none",
    table.border.bottom.style = "solid",
    table.border.bottom.color = "black",
    table.border.bottom.width = px(3),
    column_labels.border.top.color = "white",
    column_labels.border.top.width = px(3),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(3),
    data_row.padding = px(10),
    table.width=pct(100),
    container.overflow.x = F,
    container.width = 800
  ) %>% 
  tab_style(
    style = list(
      cell_text(
        size = px(38),
        weight = "normal",
        align = "left",
        color = "blue",
        font = "Bloomsbury"
      )
    ),
    locations = list(
      cells_title(groups = "title")
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(
        size = px(18),
        align = "left",
        color = "red"
      )
    ),
    locations = list(
      cells_title(groups = "subtitle")
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(
        size = px(18)
      ),
      cell_borders(
        sides = c("bottom", "top"),
        color = "green",
        weight = px(1)
      )
    ),
    locations = list(
      cells_body(gt::everything())
    )
  ) %>% 
  tab_style(
    style = list( 
      cell_text(
        size = px(14),
        color = "yellow",
        font = "Bloomsbury"
      )
    ),
    locations = list(
      cells_column_labels(everything())
    )
  ) %>% 
  tab_style(
    style = list( 
      cell_text(
        size = px(16),
        color = "red",
        font = "Bloomsbury"
      ),
      cell_borders(
        sides = c("bottom"),
        style = "solid",
        color = "darkblue",
        weight = px(2)
      )
    ),
    locations = list(
      cells_row_groups(gt::everything())
    )
  ) %>% 
  tab_style(
    style = list( 
      cell_text(
        size = px(16),
        color = "orange",
        font = "Bloomsbury"
      ),
      cell_borders(
        sides = c("bottom", "right"),
        style = "solid",
        color = "lightgrey",
        weight = px(1)
      )
    ),
    locations = list(
      cells_stub(gt::everything()),
      cells_stubhead()
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(
        font = "Bloomsbury", size = px(16), 
        color = "darkgrey")
    ),
    location = list(
      cells_body(columns = c(ID))
    )
  ) %>%
  cols_width(
    ID~pct(5),
    shelf.life~pct(17),
    cooking.time~pct(12),
    no..Of.ingredients~pct(12),
    additional.labelling~pct(30), 
    Manufacturer~pct(12)
    )%>%
  cols_align(align = "center", columns = 2:6
             )%>%
  row_group_order(groups = c("ppft", "ppmt", "spft", "spmt"))



gtsave(tab, "./06-outputfiles/table1.rtf")

library(flextable)
set_flextable_defaults(
  font.family = "Times", font.size = 10, 
  border.color = "gray")

flextable(tabledata) %>% 
  bold(part = "header") %>% 
  add_footer_lines("The 'cars' dataset")
%>% 
  save_as_docx(path = "mytable.docx")
