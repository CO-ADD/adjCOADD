$(document).ready(function () {
  $("#create_stock").on("click", function () {
    let my_modal = $("#createStockModal");

    my_modal.load("/createStock/", function (response) {
      if (response == "Permission Not Granted") {
        alert("permission!");
      } else {
        my_modal.modal("show"); // Open Modal
      }
    });
  });
});
