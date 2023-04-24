function loadModal(modalId, triggerElementSelector, url) {
    let my_modal = $(modalId);
    $(triggerElementSelector).click(function () {
     
      my_modal.load(url, function (response) {
        if (response == "Permission Not Granted") {
          alert("permission!");
        } else {
          my_modal.modal("show");
        }
      });
    });
  }