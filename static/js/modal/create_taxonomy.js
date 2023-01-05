$(document).ready(function () {
  let my_modal = $("#createTaxoModal");

  $("#createTaxo").click(function () {
    my_modal.load("/createTaxo/", function (response) {
      if (response == "Permission Not Granted") {
        alert("permission!");
      } else {
        my_modal.modal("show"); // Open Modal
      }
    });
  });
});
