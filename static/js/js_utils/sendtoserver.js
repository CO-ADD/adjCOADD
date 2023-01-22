const csrftoken = getCookie("csrftoken");
const sendToServer = (data, url) => {
  console.log("send to server");
  $.ajax({
    url: url,
    type: "POST",
    headers: { "X-CSRFToken": csrftoken },
    data: data,
  })
    .done((response) => {
      data = response;
      return data;
    })
    .fail(() => {
      console.log("Error occured");
    });
};
