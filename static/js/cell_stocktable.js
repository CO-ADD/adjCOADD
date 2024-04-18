function cellstocktable(permission, table) {

  $("#datatable_cell_batch tbody").on("click", "a.dt-control", function (e) {

    var tr = $(this).closest("tr");
    var row = table.row(tr);
    var Batch_id_a = $(this).data("value");
    var stock_create_id = "create_stock_for" + $(this).data("value").toString();
    var stock_create_modal = "create_stock_modal_for" + $(this).data("value").toString();
    if (row.child.isShown()) {
      // This row is already open - close it
      destroyCellStocktable(row);
      tr.removeClass("shown");
    } else {
      // Open this row
      // getStockData(row.data()[1]);
      createCellStocktable(row, Batch_id_a, stock_create_id, stock_create_modal, permission);
      tr.addClass("shown");
    }
  });

}