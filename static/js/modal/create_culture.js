$(document).ready(function () {
  let my_modal = $("#createCultureModal");

  $("#createCulture").click(function () {
    my_modal.load("/createCulture/", function (response) {
      if (response == "Permission Not Granted") {
        alert("permission!");
      } else {
        my_modal.modal("show"); // Open Modal
      }
    });
  });
});
