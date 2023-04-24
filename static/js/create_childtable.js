function createChild(row, Batch_id_a, stock_create_id, stock_create_modal, permission) {
// Define
const csrftoken = getCookie("csrftoken");
var defaultContent = (data, type, row) => {
  
  return permission == 'Admin' ?
    (`<button id="update_stock" class="text-success border-0"> <i class="fa-regular fa-pen-to-square"></i></button>` +
      '<div id="updateStockModal" class="modal fade" role="dialog"></div>' +
      "&nbsp; " +
      `<button id="delete_stock" class="border-0"> <i class="fa-regular fa-trash-can"></i></button>` +
      '<div id="deleteStockModal" class="modal " role="dialog"></div>')
    :
    (`<button class="text-secondary border-0" disabled> <i class="fa-regular fa-pen-to-square"></i></button>` +
      '<div id="updateStockModal" class="modal fade" role="dialog"></div>' +
      "&nbsp; " +
      `<button class="border-0" disabled> <i class="fa-regular fa-trash-can text-secondary"></i></button>` +
      '<div id="deleteStockModal" class="modal" role="dialog"></div>')

}
    // This is the table we'll convert into a DataTable
    var Batch_id = Batch_id_a//row.data()[1];
    var table = $(`<table id="StockTable${Batch_id_a}" class="display childtable" width="100%"/>`);
    // Display it the child row
    row.child(table).show();
    
    // Initialise as a DataTable
    var usersTable = table.DataTable({
      dom: "Bfrtip",
      pageLength: 1,
      ajax: {
        url: `/dorganism/stocklist/${Batch_id}`,
        type: "GET",
        headers: { "X-CSRFToken": csrftoken},
        data: { Batch_id },
      },

      columns: [
        // { title: "Stock ID", data: "stock_id" },
        {

          "title": `<a class="card-link"  id="${stock_create_id}"  style="z-index:10"><i class="bi bi-file-earmark-plus text-success"></i></a>
        <div id="${stock_create_modal}" class="modal fade" role="dialog"></div> New `,
          "data": null,
          "defaultContent": defaultContent()
        },
        { title: "Type", data: "stock_type" },
        { title: "Freezer", data: "location_freezer" },
        { title: "Rack", data: "location_rack" },
        { title: "Col", data: "location_col" },
        { title: "Slot", data: "location_slot" },

        { title: "Stock Date", data: "stock_date" },
        {
          title: "Left",
          data: "n_left",
          className: "n_left_column",//"editablerow_stock n_left_column",

        },
        { title: "Created", data: "n_created" },
        { title: "Note", data: "stock_notes" },
        { title: "Biologist", data: "biologist" },
      ],
      columnDefs: [
        {
          targets: '_all',
          createdCell: function (cell, cellData, rowData, rowIndex, colIndex) {
            if (colIndex && colIndex === 7){

              var n_left = rowData.n_left
              var trigger_modal_nleft=n_left>0? `<a class="card-link" type="button" style="float:right" data-bs-toggle="modal"
              data-bs-target="#m_${rowData.stock_id}_nleft_modal" id="m_${rowData.stock_id}_nleft">
              <i class="fa-regular fa-square-minus"></i>
              </a>` : `<button class="card-link" type="button" style="float:right" data-bs-toggle="tooltip" data-bs-placement="right"
              title="Stock left low" disabled>
              <i class="fa-regular fa-square-minus text-secondary"></i>
              </button>`
              const html =
              `
              ${trigger_modal_nleft}
                <div id="m_${rowData.stock_id}_nleft_modal" class="modal fade" tabindex="-1" aria-labelledby="exampleModalLabel"
                aria-hidden="true">
                <div class="modal-dialog">
                  <div class="modal-content">
                     
                    <div class="modal-body fw-bold text-center">
                      Click Save, n_Left will be reduced by 1

                      </div>
                      <div class="modal-footer d-flex justify-content-center align-center">
                      <button id="btn_save_${rowData.stock_id}_nleft" class="p-0 border border-0 m-auto bg-transparent text-center"><i class="fa-regular fa-floppy-disk text-success m-auto"></i></button>
                      <button type="button" class="btn-close m-auto border-0" data-bs-dismiss="modal"></button>
                      </div>
                      </div>
                   </div>
                   </div>
                   `
                   if (colIndex === 7) {
                     $(cell).attr('data-type', 'n_left');
                     if (rowData.n_left < 2) {
                $(cell).html(`left low!`)
              } else {
                $(cell).html(`${n_left}` + html)
              }
            
            }
            
            $('body').on("click", `#btn_save_${rowData.stock_id}_nleft`, function () {
              
              $.ajax({
                url: `/dorganism/updateStock/${rowData.stock_id}/`,
                type: "POST",
                headers: { "X-CSRFToken": csrftoken },
                data: { 'value': n_left.toString() },
              })
                .done((response) => {
                  n_left = response["result"]
                  console.log(n_left)
                  $(`#m_${rowData.stock_id}_nleft_modal`).modal('hide'); // Hide the modal

              // Update the cell content
                
                  cell.innerHTML = `${n_left}` + html;            
                })
                .fail((xhr, textStatus, errorThrown) => {
                  console.log("Error occured");
                  cell.innerHTML=`${xhr}, ${textStatus}, ${errorThrown}`
                });
            })

          }

          }
        }
      ],
      // select: true,
      searching: false,
      paging: false,
      info: false
    });

    $(`#StockTable${Batch_id_a} tbody`).on("click", `#delete_stock`, function () {
      var data = usersTable.row($(this).parents("tr")).data();
      let my_modal = $("#deleteStockModal");
      console.log(data);
      var pk = data["stock_id"].toString();
      console.log(pk);
      my_modal.load(`/dorganism/deleteStock/${pk}/`, function (response) {
        if (response == "Permission Not Granted") {
          alert("permission!");
        } else {
          my_modal.modal("show"); // Open Modal
        }
      });
    });

    let url_stock= `/dorganism/createStock/${Batch_id_a}/`
    
    loadModal("#" + stock_create_modal, "#" + stock_create_id, url_stock);
    $(`#StockTable${Batch_id_a} tbody`).on("click", `#update_stock`, function () {
      var data = usersTable.row($(this).parents("tr")).data();
      let my_modal = $("#updateStockModal");
      console.log(data);
      var pk = data["stock_id"].toString();
      console.log(pk);
      my_modal.load(`/dorganism/updateStock/${pk}/`, function (response) {
        if (response == "Permission Not Granted") {
          alert("permission!");
        } else {
          my_modal.modal("show"); // Open Modal
        }
      });
    });
  }