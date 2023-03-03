
// $(document).ready(function () {


// });
(function ($, id_modal, id_trigger, id_drug) {
    let my_modal = $(id_modal);

    $(id_trigger).click(function () {
        console.log("click structrue")
        my_modal.load(`/drug_detail_structure/${id_drug}`, function (response) {
            if (response == "Permission Not Granted") {
                alert("permission!");
            } else {
                my_modal.modal("show"); // Open Modal
            }
        });
    });

}(jQuery));