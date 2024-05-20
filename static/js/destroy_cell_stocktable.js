function destroyCellStocktable(row, Batch_id_a) {
    var table = $("table", row.child());
    table.detach();
    table.DataTable().destroy();

    // And then hide the row
    row.child.remove();
  }
