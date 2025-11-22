/**
 * Notification Utility - Unified toast notification system
 * 统一的消息通知系统
 */
import React from 'react'
import { message, notification } from 'antd'
import type { ArgsProps as MessageArgsProps } from 'antd/es/message'
import type { ArgsProps as NotificationArgsProps } from 'antd/es/notification'
import {
  CheckCircleOutlined,
  CloseCircleOutlined,
  ExclamationCircleOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons'

// Configure global message settings
message.config({
  top: 80,
  duration: 3,
  maxCount: 3,
  rtl: false,
})

// Configure global notification settings
notification.config({
  placement: 'topRight',
  top: 80,
  duration: 4.5,
  rtl: false,
})

/**
 * Toast notifications - Simple short messages
 */
export const toast = {
  success: (content: string, duration?: number) => {
    message.success({
      content,
      duration,
      icon: <CheckCircleOutlined style={{ color: '#10b981' }} />,
    })
  },

  error: (content: string, duration?: number) => {
    message.error({
      content,
      duration,
      icon: <CloseCircleOutlined style={{ color: '#ef4444' }} />,
    })
  },

  warning: (content: string, duration?: number) => {
    message.warning({
      content,
      duration,
      icon: <ExclamationCircleOutlined style={{ color: '#f59e0b' }} />,
    })
  },

  info: (content: string, duration?: number) => {
    message.info({
      content,
      duration,
      icon: <InfoCircleOutlined style={{ color: '#2563eb' }} />,
    })
  },

  loading: (content: string) => {
    return message.loading(content, 0)
  },
}

/**
 * Notifications - Detailed messages with title and description
 */
export const notify = {
  success: (title: string, description?: string, options?: NotificationArgsProps) => {
    notification.success({
      message: title,
      description,
      icon: <CheckCircleOutlined style={{ color: '#10b981' }} />,
      ...options,
    })
  },

  error: (title: string, description?: string, options?: NotificationArgsProps) => {
    notification.error({
      message: title,
      description,
      icon: <CloseCircleOutlined style={{ color: '#ef4444' }} />,
      ...options,
    })
  },

  warning: (title: string, description?: string, options?: NotificationArgsProps) => {
    notification.warning({
      message: title,
      description,
      icon: <ExclamationCircleOutlined style={{ color: '#f59e0b' }} />,
      ...options,
    })
  },

  info: (title: string, description?: string, options?: NotificationArgsProps) => {
    notification.info({
      message: title,
      description,
      icon: <InfoCircleOutlined style={{ color: '#2563eb' }} />,
      ...options,
    })
  },
}

/**
 * Specialized notifications for common scenarios
 */
export const notifications = {
  // Success scenarios
  saveSuccess: () => toast.success('保存成功'),
  updateSuccess: () => toast.success('更新成功'),
  deleteSuccess: () => toast.success('删除成功'),
  createSuccess: (itemName: string) => toast.success(`${itemName}创建成功`),
  uploadSuccess: () => toast.success('上传成功'),
  submitSuccess: () => toast.success('提交成功'),

  // Error scenarios
  saveError: () => toast.error('保存失败，请重试'),
  updateError: () => toast.error('更新失败，请重试'),
  deleteError: () => toast.error('删除失败，请重试'),
  createError: (itemName: string) => toast.error(`${itemName}创建失败，请重试`),
  uploadError: () => toast.error('上传失败，请重试'),
  submitError: () => toast.error('提交失败，请重试'),
  networkError: () => toast.error('网络连接失败，请检查网络'),
  serverError: () => toast.error('服务器错误，请稍后重试'),
  unauthorized: () => toast.error('您没有权限执行此操作'),
  notFound: () => toast.error('请求的资源不存在'),
  validationError: (field: string) => toast.error(`${field}验证失败`),

  // Loading scenarios
  saving: () => toast.loading('保存中...'),
  loading: () => toast.loading('加载中...'),
  uploading: () => toast.loading('上传中...'),
  processing: () => toast.loading('处理中...'),

  // Task execution notifications
  taskStarted: (taskName: string) =>
    notify.info('任务已启动', `${taskName}已开始执行，您可以在任务列表中查看进度`),
  taskCompleted: (taskName: string) =>
    notify.success('任务完成', `${taskName}已成功完成`),
  taskFailed: (taskName: string, error?: string) =>
    notify.error('任务失败', error || `${taskName}执行失败，请查看日志`),

  // Pipeline execution notifications
  pipelineStarted: (pipelineName: string) =>
    notify.info(
      '管道已启动',
      `${pipelineName}已开始运行，预计需要几分钟到几小时`,
      { duration: 6 }
    ),
  pipelineCompleted: (pipelineName: string) =>
    notify.success('管道完成', `${pipelineName}已成功完成，可以查看结果`),
  pipelineFailed: (pipelineName: string) =>
    notify.error('管道失败', `${pipelineName}运行失败，请检查参数和输入文件`),

  // File management notifications
  fileUploaded: (filename: string) =>
    toast.success(`文件 "${filename}" 上传成功`),
  fileDeleted: (filename: string) =>
    toast.success(`文件 "${filename}" 已删除`),
  fileDownloaded: (filename: string) =>
    toast.success(`文件 "${filename}" 开始下载`),

  // Form validation
  fillRequired: () => toast.warning('请填写所有必填项'),
  invalidFormat: (field: string) => toast.warning(`${field}格式不正确`),
  confirmDelete: (itemName: string) =>
    notify.warning(
      '确认删除',
      `您确定要删除 "${itemName}" 吗？此操作不可撤销`,
      { duration: 0 }
    ),
}

/**
 * API error handler - Convert API errors to user-friendly messages
 */
export const handleApiError = (error: any, defaultMessage = '操作失败') => {
  if (!error) {
    toast.error(defaultMessage)
    return
  }

  // Network error
  if (!error.response) {
    notifications.networkError()
    return
  }

  // HTTP status codes
  const status = error.response?.status
  const message = error.response?.data?.message || error.message

  switch (status) {
    case 400:
      toast.error(message || '请求参数错误')
      break
    case 401:
      toast.error('您的登录已过期，请重新登录')
      // Redirect to login page
      setTimeout(() => {
        window.location.href = '/login'
      }, 1500)
      break
    case 403:
      notifications.unauthorized()
      break
    case 404:
      notifications.notFound()
      break
    case 422:
      toast.error(message || '数据验证失败')
      break
    case 500:
      notifications.serverError()
      break
    case 503:
      toast.error('服务暂时不可用，请稍后重试')
      break
    default:
      toast.error(message || defaultMessage)
  }
}

export default {
  toast,
  notify,
  notifications,
  handleApiError,
}
