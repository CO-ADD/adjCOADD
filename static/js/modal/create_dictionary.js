$(document).ready(function () {
  let my_modal = $("#createDictModal");

  $("#createDict").click(function () {
    my_modal.load("/dict_create/", function (response) {
      if (response == "Permission Not Granted") {
        alert("permission!");
      } else {
        my_modal.modal("show"); // Open Modal
      }
    });
  });
});
