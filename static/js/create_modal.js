const loadModal = (modalId, triggerElementSelector, url, csrftoken, formSelector = null) => {
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

  // Add the form submission handler here, only when formSelector is provided
  if (formSelector) {
    $(document).on("submit", formSelector, function (e) {
      e.preventDefault(); // Prevent the default form submission

      // Serialize the form data
      const formData = $(this).serialize();

      // Make an AJAX POST request to your view
      $.ajax({
        url: url, // Use the same URL to submit the form
        type: "POST",
        data: formData,
        headers: { "X-CSRFToken": csrftoken },
        success: function (response) {
          // If successful, close the modal and refresh the page or update the content.
          my_modal.modal("hide");
          location.reload();
        },
        error: function (xhr) {
          // If there is an error, update the form with the errors
          my_modal.modal("show");
          const errors = xhr.responseJSON.errors;
          for (const field in errors) {
            const errorMessage = errors[field];
            $(`#id_${field}`)
              .siblings(".error-message")
              .text(errorMessage);
          }
        },
      });
    });
  }
};
