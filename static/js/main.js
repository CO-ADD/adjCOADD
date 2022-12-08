// custom javascript

$(document).ready(() => {
  console.log('Sanity Check!');
});

$('.button').on('click', function () {
  console.log("button clicked")
  const filepath = $("#filepath").text() ? $("#filepath").text() : "none"
  console.log(filepath)
  $.ajax({
    url: '/tasks/',
    data: { type: $(this).data('type'), filepath: filepath.toString() },
    method: 'POST',
  })
    .done((res) => {
      console.log(res);
      // const html = `
      // <tr>
      //   <td>${res.task_id}</td>
      //   <td>${res.task_status}</td>
      //   <td>${res.task_result}</td>
      // </tr>`
      // $('#tasks').prepend(html);
      getStatus(res);
    })
    .fail((err) => {
      console.log(err);
    });
});

function getStatus(res) {
  console.log('this')
  console.log(res)
  $.ajax({
    url: `/tasks/${res.task_id}/`,
    method: 'GET'
  })
    .done((res) => {
      console.log(`getstatus ${res.task_status}`)
      const html = `
      <tr>
        <td>${res.task_id}</td>
        <td>${res.task_status}</td>
        <td>${res.task_result}</td>
      </tr>`
      $('#tasks').prepend(html);

      const taskStatus = res.task_status;

      if (taskStatus === 'pass' || taskStatus === 'warning') {
        save_object(res.task_id, file_url)
      } else {
        console.log('error List')
      }
      // setTimeout(function () {
      //   getStatus(res.task_id);
      // }, 1000);
    })
    .fail((err) => {
      console.log(err)
    });



}

function save_object(filepath) {
  $.ajax({
    url: '/tasks/',
    data: { type: $(this).data('type') },
    method: 'POST',
  })
    .done((res) => {
      console.log(res)
      getStatus(res);
    })
    .fail((err) => {
      console.log(err);
    });

}
