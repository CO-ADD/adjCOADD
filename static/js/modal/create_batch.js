$(document).ready(function () {
  let my_modal = $("#createBatchModal");

  $("#createBatch").click(function () {
    my_modal.load("/createBatch/", function (response) {
      if (response == "Permission Not Granted") {
        alert("permission!");
      } else {
        my_modal.modal("show"); // Open Modal
      }
    });
  });
});
