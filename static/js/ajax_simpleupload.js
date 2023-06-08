function ajaxSimpleUpload(inputSelector,inputfield, url, fileNameInputSelector, fileTypeInputSelector) {

    const csrftoken = getCookie("csrftoken");
    $(inputSelector).on('change', function () {
        console.log(url)
        var fileInput = $(inputSelector)[0];
        var file = fileInput.files[0];
        var formData = new FormData();
        formData.append(inputfield, file);
        $.ajax({
            url: url,
            type: 'POST',
            data: formData,
            processData: false,
            contentType: false,
            headers: { "X-CSRFToken": csrftoken},
            success: function (data) {
                $(fileNameInputSelector).val(data.name);
                $(fileTypeInputSelector).val(data.type);
            },
        });
    });
}